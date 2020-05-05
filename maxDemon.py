#!/usr/bin/env python
import os
from kivy.app import App
from kivy.uix.widget import Widget
from kivy.properties import ObjectProperty

class MAXDEMONApp(App):
#    kv_directory = 'kivy_templates'
    def build(self):
        return MyLayout()
    global s
    s = 'off'
    
class MyLayout(Widget):
    
    fileML = open("MLmethods.txt","w")
    fileML.write("MLmethods\n")
    fileML.write("bnp\n")
    fileML.write("dist\n")
    fileML.write("kern\n")
    fileML.write("ens\n")
    fileML.write("shape\n")
    fileML.close()
    
    
    
    # define buttons and actions
     
    def ml1(self):
            print("KNN ML deselected")
            fileML = open("MLmethods.txt","a+")
            fileML.writelines("no_bnp\n")
            fileML.close()
    def ml2(self):
            print("probablistic ML deselected")
            fileML = open("MLmethods.txt","a+")
            fileML.writelines("no_dist\n")
            fileML.close()
    def ml3(self):
            print("SVM ML deselected")
            fileML = open("MLmethods.txt","a+")
            fileML.writelines("no_kern\n")
            fileML.close()
    def ml4(self):
            print("ensemble ML deselected")
            fileML = open("MLmethods.txt","a+")
            fileML.writelines("no_ens\n")
            fileML.close()
    def ml5(self):
            print("shape analysis selected")
            fileML = open("MLmethods.txt","a+")
            fileML.writelines("no_shape\n")
            fileML.close()
            global s
            s = 'on'
    def btn1(self):
            print ("Please enter type of structure (1=protein | 2=DNA+protein | 3=protein+ligand)")
            structure = input("Enter 1,2 or 3: ")
            print("running validation and variant MD simulations") 
            if structure == 1:
                cmd = 'perl GUI_MLMD_DROIDS.pl'
                os.system(cmd)
            if structure == 2:
                cmd = 'perl GUI_MLMD_DROIDSdp.pl'
                os.system(cmd)
            if structure == 3:
                cmd = 'perl GUI_MLMD_DROIDSlp.pl'
                os.system(cmd)
    def btn2(self):
            print("making control file and processing training data sets") 
            cmd = 'perl GUI_ML_DROIDSbtn2.pl'
            os.system(cmd)
    def btn3(self):
            if s == 'on':
                print("running stacked machine learning model") 
                cmd = 'perl GUI_ML_DROIDSbtn3.pl'
                os.system(cmd)
            if s == 'off':
                print("running stacked machine learning model") 
                cmd = 'perl GUI_ML_DROIDSbtn3_fluxOnly.pl'
                os.system(cmd)    
    def btn4(self):
            print("running canonical correlation and variant impact analyses") 
            cmd = 'perl GUI_ML_DROIDSbtn4.pl'
            os.system(cmd)
    def btn5(self):
            print("showing conserved dynamics on protein structure") 
            cmd = 'perl GUI_ML_DROIDSbtn5.pl'
            os.system(cmd)        
    def btn6(self):
            print("rendering movies") 
            cmd = 'perl GUI_ML_DROIDSbtn6.pl'
            os.system(cmd)        
    def btn7(self):
            print("showing variant impacts on protein structure") 
            cmd = 'perl GUI_ML_DROIDSbtn7.pl'
            os.system(cmd)        
    def btn8(self):
            print("playing movies of learners on protein structure") 
            cmd = 'perl GUI_ML_DROIDSbtn8.pl'
            os.system(cmd)        
    
if __name__ == '__main__':
    MAXDEMONApp().run()
