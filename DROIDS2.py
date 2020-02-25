#!/usr/bin/env python
import os
from kivy.app import App
from kivy.uix.widget import Widget
from kivy.properties import ObjectProperty

class DROIDS2App(App):
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
            print("running DROIDS on single GPU - compare mutations on single protein") 
            cmd = 'perl GUI_START_DROIDSsdm.pl'
            os.system(cmd)
        if self.gpu == 2:
            print("running DROIDS on dual GPU - compare mutations on single protein") 
            cmd = 'perl GUI_START_DROIDSsdm_dualGPU.pl'
            os.system(cmd)
    def btn2(self):
        if self.gpu == 1:
            print("running DROIDS on single GPU - compare mutations on DNA-protein interaction") 
            cmd = 'perl GUI_START_DROIDSdp2.pl'
            os.system(cmd)
        if self.gpu == 2:
            print("running DROIDS on dual GPU - compare mutations on DNA-protein interaction") 
            cmd = 'perl GUI_START_DROIDSdp2_dualGPU.pl'
            os.system(cmd)
    def btn3(self):
        if self.gpu == 1:
            print("running DROIDS on single GPU - compare mutations on protein-ligand interaction") 
            cmd = 'perl GUI_START_DROIDSlp2.pl'
            os.system(cmd)
        if self.gpu == 2:
            print("running DROIDS on dual GPU - compare mutations on protein-ligand interaction") 
            cmd = 'perl GUI_START_DROIDSlp2_dualGPU.pl'
            os.system(cmd)
    


if __name__ == '__main__':
    DROIDS2App().run()
