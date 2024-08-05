import numpy as np
import xarray as xr


# from Stiperski and Calaf 2018

def triangles(state, perc):
    '''From the oordinates of the 6 triangles with vertex at the barycenter of the equilateral triangle of side=1, returns the two orresponding to the asymptotic state including the percentage cut'''
    #triangles are row = vertex number, col = x,y
    bary = np.array([0.5,np.sqrt(3)/6])
    A = np.array([0,0])
    B = np.array([0.5*np.cos(np.pi/180*60), 0.5*np.sin(np.pi/180*60)])
    C = np.array([0.5, np.sqrt(3)/2])
    D = np.array([1-B[0],B[1]])
    E = np.array([1,0])
    F = np.array([0.5,0])
    
    if state=='1c':
        bary2 = [1-bary[0]*perc, bary[1]*perc]
        triangle1 = [E, [1-B[0]*perc,D[1]*perc], bary2]
        triangle2 = [E, F*(2-perc), bary2]
    
    elif state=='2c':
        bary2 = bary*perc
        triangle1 = [A, B*perc, bary2]
        triangle2 = [A, F*perc, bary2]
        
    elif state=='3c':
        bary2 = [bary[0],C[1]-(C[1]-bary[1])*perc]
        triangle1 = [C, B*(2-perc), bary2]
        triangle2 = [C, [1-B[0]*(2-perc),D[1]*(2-perc)], bary2]
    
    return np.array(triangle1), np.array(triangle2)

def isin_triangle(triangle, x, y):
    
    x1,y1 = triangle[0]
    x2,y2 = triangle[1]
    x3,y3 = triangle[2]
    
    c1 = (x2-x1)*(y-y1)-(y2-y1)*(x-x1)
    c2 = (x3-x2)*(y-y2)-(y3-y2)*(x-x2)
    c3 = (x1-x3)*(y-y3)-(y1-y3)*(x-x3)
    
    return ((c1<0) & (c2<0) & (c3<0)) | ((c1>0) & (c2>0) & (c3>0))
    
def select_asymptotic(ds, state, perc=0.7):

    triangle1, triangle2 = triangles(state, perc)
    x = ds.xb
    y = ds.yb
    belongs = isin_triangle(triangle1,x,y) | isin_triangle(triangle2,x,y)

    return ds.where(belongs)