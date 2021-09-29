# -*- coding: utf-8 -*-
"""
Created on Fri Sep 24 10:33:20 2021

@author: jks68 and sunny who am I kidding
"""

from os import listdir
from os.path import isfile, join
import os



def sunnyfunction(directory,searchterm):
    allfiles = []
    for path, subdirs, files in os.walk(directory):
        for name in files:
            allfiles.append(str(os.path.join(path, name)))
   
    onlypython = [f for f in allfiles if f[-3:] == '.py']
    print(str(len(onlypython)) + ' files found')
   
    a =0
   
    for i in range(len(onlypython)):
        f = open(onlypython[i], "r")
        code = f.read()
       
        if searchterm in code:
            print('match found: ' + str(onlypython[i]))
            a=1
    if a ==0:
        print('no match found')
   


if __name__ == '__main__':
       
    directory = "C:\Users\jks68\Documents\GitHub"
    searchterm = 'AndorControlUI' #put the part of a code that you wanna find here!!!!!!!!!!!!!!!!!!
    sunnyfunction(directory,searchterm)