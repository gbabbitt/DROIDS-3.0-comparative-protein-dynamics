#!/usr/bin/env python
import os
from kivy.app import App
from kivy.uix.widget import Widget
from kivy.properties import ObjectProperty

print("Welcome to DROIDS- Detecting Relative Outlier Impacts in Dynamic Simulations- analytical engine and visual toolbox for functional evolutionary comparison of molecular dynamic simulation")
cmd = 'gedit READMEv3.0.md'
os.system(cmd)
print("finding paths for paths.ctl") 
cmd = 'perl PATHS.pl'
os.system(cmd)


class DROIDSApp(App):
#    kv_directory = 'kivy_templates'
    def build(self):
        return MyLayout()
   
class MyLayout(Widget):
    
      
    # define buttons and actions
    def btn1(self):
        print("running DROIDS - direct comparative analysis") 
        cmd = 'python DROIDS1.py'
        os.system(cmd)
    def btn2(self):
        print("running DROIDS - mutant model comparison") 
        cmd = 'python DROIDS2.py'
        os.system(cmd)
    def btn3(self):
        print("running DROIDS + maxDemon - functional variant analysis") 
        cmd = 'python DROIDS3.py'
        os.system(cmd)
       


if __name__ == '__main__':
    DROIDSApp().run()
