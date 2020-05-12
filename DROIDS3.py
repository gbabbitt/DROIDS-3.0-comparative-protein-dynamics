#!/usr/bin/env python
import os
from kivy.app import App
from kivy.uix.widget import Widget
from kivy.properties import ObjectProperty

class DROIDS3App(App):
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
            print("running DROIDS on single GPU - single protein analysis") 
            cmd = 'perl GUI_START_DROIDSss1.pl'
            os.system(cmd)
        if self.gpu == 2:
            print("running DROIDS on dual GPU - single protein analysis") 
            cmd = 'perl GUI_START_DROIDSss1_dualGPU.pl'
            os.system(cmd)
    def btn2(self):
        if self.gpu == 1:
            print("running DROIDS on single GPU - protein-DNA analysis") 
            cmd = 'perl GUI_START_DROIDSdp1.pl'
            os.system(cmd)
        if self.gpu == 2:
            print("running DROIDS on dual GPU - protein-DNA analysis") 
            cmd = 'perl GUI_START_DROIDSdp1_dualGPU.pl'
            os.system(cmd)
    def btn3(self):
        if self.gpu == 1:
            print("running DROIDS on single GPU - protein-ligand analysis") 
            cmd = 'perl GUI_START_DROIDSlp1.pl'
            os.system(cmd)
        if self.gpu == 2:
            print("running DROIDS on dual GPU - protein-ligand analysis") 
            cmd = 'perl GUI_START_DROIDSlp1_dualGPU.pl'
            os.system(cmd)
    def btn4(self):
        if self.gpu == 1:
            print("running DROIDS on single GPU - protein-protein analysis") 
            cmd = 'perl GUI_START_DROIDSpp.pl'
            os.system(cmd)
        if self.gpu == 2:
            print("running DROIDS on dual GPU - protein-protein analysis") 
            cmd = 'perl GUI_START_DROIDSpp_dualGPU.pl'
            os.system(cmd)


if __name__ == '__main__':
    DROIDS3App().run()
