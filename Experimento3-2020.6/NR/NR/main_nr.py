from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtGui import QFont
from PyQt5.QtWidgets import *
import sys
from GUI.myui import * #importando somendo o arquivo .py da pasta GUI
from calc_nr import calcTputNR  #importando o arquivo com a função e dicionários



class calcMaxTput(QtWidgets.QMainWindow , Ui_MainWindow):

    def __init__(self, parent = None):
        super(calcMaxTput, self).__init__(parent)
        self.setupUi(self)
        self.setWindowTitle('Calculadora 5G NR R15')
        self.pushButton.clicked.connect(self.button_action)


    def button_action(self):

        mimo = int(self.comboBox_mimo.currentText())
        carrierAgregation  = int(self.comboBox_carrierAgreg.currentText())
        modulationOrder = int(self.comboBox_morOrder.currentText())
        numerology = int(self.comboBox_Numerology.currentText())
        overhead = float(self.comboBox_overhead.currentText())
        numRBs = int(self.prbs.text())
        maxTput = calcTputNR(carrierAgregation, modulationOrder , mimo, numerology, numRBs , overhead)
        # maxTputM = f'{maxTput:.5f}'
        # maxTputG = f'{maxTput/1e3:.5f}'
        # self.maxTputMhz.display(maxTputM)
        # self.maxTputGHZ.display(maxTputG)
        self.maxTputMhz.setAlignment(QtCore.Qt.AlignCenter)
        self.maxTputMhz.setFont(QFont('Arial',14))
        self.maxTputMhz.setText(f'{maxTput:.5f} MHz') 
        self.maxTputGhz.setAlignment(QtCore.Qt.AlignCenter)
        self.maxTputGhz.setFont(QFont('Arial',14))
        self.maxTputGhz.setText(f'{maxTput/1e3:.5f} GHz') 


if __name__ == '__main__':
    qt = QApplication(sys.argv)
    calcmaxtput = calcMaxTput()
    calcmaxtput.show()
    qt.exec_()


    

