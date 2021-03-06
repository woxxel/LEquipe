import numpy as np

def func_const(x,m):
  return m

def func_linear_zero(x,m):
  return m*x

def func_linear(x,m,b):
  return m*x + b

def func_exp(x,a,b):
  return b*np.exp(a*x)

def func_double_exp(x,a1,a2,b1,b2):
  return b1*np.exp(a1*x)+b2*np.exp(a2*x)

def func_exp_exp(x,a,b):
  return np.exp(b*np.exp(a*x))

def func_exp_intersect(x,a1,a2,b1,b2):
  output = np.zeros(len(x))
  
  intersect = np.log(b2/b1)/(a1-a2)
  #print "inter: ", intersect
  idx_intersect = np.argmin(abs(x-intersect))
  
  output[:idx_intersect] = b1*np.exp(a1*x[:idx_intersect])
  output[idx_intersect:] = b2*np.exp(a2*x[idx_intersect:])
  
  return output
  
def func_poly(x,a,b):
  return b*x**a
  
def func_log(x,a,b):
  return a*np.log(x) + b