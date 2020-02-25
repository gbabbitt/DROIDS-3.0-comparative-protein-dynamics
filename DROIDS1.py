#!/usr/bin/env python
import os
from kivy.app import App
from kivy.uix.widget import Widget
from kivy.properties import ObjectProperty

class DROIDS1App(App):
#    kv_directory = 'kivy_templates'
    def build(self):
        return MyLayout()
   
class MyLayout(Widget):
    
      
    # define buttons and actions
    def btn0(self):
            print("dual GPU mode selected")
            self.gpu = 2
            return self.gpu
    def btn00(self):
            print("single GPU mode selected")
            self.gpu = 1
            return self.gpu    
    def btn1(self):
        if self.gpu == 1:
            print("running DROIDS on single GPU - compare single protein at two temperatures") 
            cmd = 'perl GUI_START_DROIDSss.pl'
            os.system(cmd)
        if self.gpu == 2:
            print("running DROIDS on dual GPU - compare single protein at two temperatures") 
            cmd = 'perl GUI_START_DROIDSss_dualGPU.pl'
            os.system(cmd)
    def btn2(self):
        if self.gpu == 1:
            print("running DROIDS on single GPU - compare two evolutionary homologs") 
            cmd = 'perl GUI_START_DROIDSed.pl'
            os.system(cmd)
        if self.gpu == 2:
            print("running DROIDS on dual GPU - compare two evolutionary homologs") 
            cmd = 'perl GUI_START_DROIDSed_dualGPU.pl'
            os.system(cmd)
    def btn3(self):
        if self.gpu == 1:
            print("running DROIDS on single GPU - compare two DNA-protein binding interactions") 
            cmd = 'perl GUI_START_DROIDSdp3.pl'
            os.system(cmd)
        if self.gpu == 2:
            print("running DROIDS on dual GPU - compare two DNA-protein binding interactions") 
            cmd = 'perl GUI_START_DROIDSdp3_dualGPU.pl'
            os.system(cmd)
    def btn4(self):
        if self.gpu == 1:
            print("running DROIDS on single GPU - compare two protein-ligand interactions") 
            cmd = 'perl GUI_START_DROIDSlp3.pl'
            os.system(cmd)
        if self.gpu == 2:
            print("running DROIDS on dual GPU - compare two protein-ligand interactions") 
            cmd = 'perl GUI_START_DROIDSlp3_dualGPU.pl'
            os.system(cmd)


if __name__ == '__main__':
    DROIDS1App().run()
