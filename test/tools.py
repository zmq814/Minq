## test version
import os
from numpy import *

class Mydata:
  def __init__{self,data}:
    self.file = data
    self.basename = os.path.basename(data)
    self.path = os.path.dir(data)
  
  def addcol(self,data):
    return mean(data, axis=1)
  
