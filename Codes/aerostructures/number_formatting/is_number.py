# -*- coding: utf-8 -*-
"""
@author: © Joan Mas Colomer
"""

#Function that checks whether a string can be converted into a float
def isfloat(value):
  try:
    float(value)
    return True
  except ValueError:
    return False

#Function that checks whether a string can be converted into an integer
def isint(value):
  try:
    int(value)
    return True
  except ValueError:
    return False
