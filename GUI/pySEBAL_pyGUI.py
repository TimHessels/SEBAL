# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'pySEBAL_GUI.ui'
#
# Created by: PyQt5 UI code generator 5.6
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_pySEBALGUI(object):
    def setupUi(self, pySEBALGUI):
        pySEBALGUI.setObjectName("pySEBALGUI")
        pySEBALGUI.resize(588, 315)
        pySEBALGUI.setMinimumSize(QtCore.QSize(588, 290))
        icon = QtGui.QIcon()
        icon.addPixmap(QtGui.QPixmap(":/Logo/tree.png"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        pySEBALGUI.setWindowIcon(icon)
        self.ButtonBox_OKCancel_SEBAL = QtWidgets.QDialogButtonBox(pySEBALGUI)
        self.ButtonBox_OKCancel_SEBAL.setGeometry(QtCore.QRect(230, 250, 341, 32))
        self.ButtonBox_OKCancel_SEBAL.setOrientation(QtCore.Qt.Horizontal)
        self.ButtonBox_OKCancel_SEBAL.setStandardButtons(QtWidgets.QDialogButtonBox.Cancel|QtWidgets.QDialogButtonBox.Ok)
        self.ButtonBox_OKCancel_SEBAL.setCenterButtons(False)
        self.ButtonBox_OKCancel_SEBAL.setObjectName("ButtonBox_OKCancel_SEBAL")
        self.pushButton = QtWidgets.QPushButton(pySEBALGUI)
        self.pushButton.setGeometry(QtCore.QRect(470, 130, 93, 31))
        self.pushButton.setObjectName("pushButton")
        self.lineEdit = QtWidgets.QLineEdit(pySEBALGUI)
        self.lineEdit.setGeometry(QtCore.QRect(30, 130, 431, 31))
        self.lineEdit.setObjectName("lineEdit")
        self.lineEdit_2 = QtWidgets.QLineEdit(pySEBALGUI)
        self.lineEdit_2.setEnabled(False)
        self.lineEdit_2.setGeometry(QtCore.QRect(210, 200, 71, 22))
        self.lineEdit_2.setAutoFillBackground(False)
        self.lineEdit_2.setReadOnly(True)
        self.lineEdit_2.setObjectName("lineEdit_2")
        self.spinBox_2 = QtWidgets.QSpinBox(pySEBALGUI)
        self.spinBox_2.setGeometry(QtCore.QRect(290, 200, 42, 22))
        self.spinBox_2.setMinimum(2)
        self.spinBox_2.setMaximum(9999)
        self.spinBox_2.setObjectName("spinBox_2")
        self.label = QtWidgets.QLabel(pySEBALGUI)
        self.label.setGeometry(QtCore.QRect(320, 0, 261, 111))
        self.label.setObjectName("label")
        self.label_2 = QtWidgets.QLabel(pySEBALGUI)
        self.label_2.setGeometry(QtCore.QRect(20, 0, 291, 111))
        self.label_2.setObjectName("label_2")
        self.groupBox = QtWidgets.QGroupBox(pySEBALGUI)
        self.groupBox.setGeometry(QtCore.QRect(10, 110, 561, 61))
        self.groupBox.setObjectName("groupBox")
        self.groupBox_2 = QtWidgets.QGroupBox(pySEBALGUI)
        self.groupBox_2.setGeometry(QtCore.QRect(10, 180, 561, 61))
        self.groupBox_2.setObjectName("groupBox_2")
        self.lineEdit_3 = QtWidgets.QLineEdit(self.groupBox_2)
        self.lineEdit_3.setEnabled(False)
        self.lineEdit_3.setGeometry(QtCore.QRect(60, 20, 71, 22))
        self.lineEdit_3.setAutoFillBackground(False)
        self.lineEdit_3.setReadOnly(True)
        self.lineEdit_3.setObjectName("lineEdit_3")
        self.spinBox = QtWidgets.QSpinBox(self.groupBox_2)
        self.spinBox.setGeometry(QtCore.QRect(140, 20, 42, 22))
        self.spinBox.setMinimum(2)
        self.spinBox.setMaximum(9999)
        self.spinBox.setObjectName("spinBox")
        self.groupBox_2.raise_()
        self.groupBox.raise_()
        self.ButtonBox_OKCancel_SEBAL.raise_()
        self.pushButton.raise_()
        self.lineEdit.raise_()
        self.lineEdit_2.raise_()
        self.spinBox_2.raise_()
        self.label.raise_()
        self.label_2.raise_()

        self.retranslateUi(pySEBALGUI)
        self.ButtonBox_OKCancel_SEBAL.accepted.connect(pySEBALGUI.accept)
        self.ButtonBox_OKCancel_SEBAL.rejected.connect(pySEBALGUI.reject)
        QtCore.QMetaObject.connectSlotsByName(pySEBALGUI)

    def retranslateUi(self, pySEBALGUI):
        _translate = QtCore.QCoreApplication.translate
        pySEBALGUI.setWindowTitle(_translate("pySEBALGUI", "pySEBAL"))
        self.ButtonBox_OKCancel_SEBAL.setToolTip(_translate("pySEBALGUI", "<html><head/><body><p>Push OK to run SEBAL or Cancel to close Tab</p></body></html>"))
        self.pushButton.setText(_translate("pySEBALGUI", "Browse"))
        self.lineEdit_2.setText(_translate("pySEBALGUI", "End Row:"))
        self.label.setText(_translate("pySEBALGUI", "<html><head/><body><p><img src=\":/Logo/ihe-delft-logo-resize.jpg\"/></p></body></html>"))
        self.label_2.setText(_translate("pySEBALGUI", "<html><head/><body><p><img src=\":/Logo/logoSEBAL-resize.png\"/></p></body></html>"))
        self.groupBox.setTitle(_translate("pySEBALGUI", "Select pySEBAL input Excel file"))
        self.groupBox_2.setTitle(_translate("pySEBALGUI", "Run Excel rows:"))
        self.lineEdit_3.setText(_translate("pySEBALGUI", "Start Row:"))

import images_rc

if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    pySEBALGUI = QtWidgets.QDialog()
    ui = Ui_pySEBALGUI()
    ui.setupUi(pySEBALGUI)
    pySEBALGUI.show()
    sys.exit(app.exec_())

